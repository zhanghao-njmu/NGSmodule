import { useState } from 'react'
import { useNavigate } from 'react-router-dom'
import { Badge, Dropdown, List, Button, Typography, Space, Empty, Spin } from 'antd'
import {
  BellOutlined,
  CheckCircleOutlined,
  InfoCircleOutlined,
  WarningOutlined,
  CloseCircleOutlined,
  CheckOutlined,
  DeleteOutlined,
  EyeOutlined,
} from '@ant-design/icons'
import type { Notification } from '@/types/notification'
import {
  useNotifications,
  useUnreadNotificationCount,
  useMarkNotificationRead,
  useMarkAllNotificationsRead,
  useDeleteNotification,
} from '@/hooks/queries'
import { DesignTokens } from '@/styles/design-tokens'
import './NotificationDropdown.css'

const { Text } = Typography

interface NotificationDropdownProps {
  className?: string
}

export const NotificationDropdown: React.FC<NotificationDropdownProps> = ({ className }) => {
  const navigate = useNavigate()
  const [dropdownVisible, setDropdownVisible] = useState(false)

  // Server state via TanStack Query — automatically refetches when cache
  // is invalidated by the realtime WebSocket layer (useNotificationStream).
  const { data: listResponse, isLoading } = useNotifications({ pageSize: 5 } as any)
  const { data: unread } = useUnreadNotificationCount()

  const markRead = useMarkNotificationRead()
  const markAllRead = useMarkAllNotificationsRead()
  const deleteOne = useDeleteNotification()

  const notifications = (listResponse?.items ?? []) as Notification[]
  const unreadCount = unread?.count ?? 0

  const handleMarkAsRead = (id: string, event: React.MouseEvent) => {
    event.stopPropagation()
    markRead.mutate(id)
  }

  const handleMarkAllAsRead = () => {
    markAllRead.mutate()
  }

  const handleDelete = (id: string, event: React.MouseEvent) => {
    event.stopPropagation()
    deleteOne.mutate(id)
  }

  const handleNotificationClick = (notification: Notification) => {
    if (!notification.read) {
      markRead.mutate(notification.id)
    }
    if ((notification as any).action_url || notification.actionUrl) {
      navigate(((notification as any).action_url || notification.actionUrl) as string)
      setDropdownVisible(false)
    }
  }

  const getNotificationIcon = (type: string) => {
    switch (type) {
      case 'success':
        return <CheckCircleOutlined style={{ color: DesignTokens.colors.success.main }} />
      case 'error':
        return <CloseCircleOutlined style={{ color: DesignTokens.colors.error.main }} />
      case 'warning':
        return <WarningOutlined style={{ color: DesignTokens.colors.warning.main }} />
      default:
        return <InfoCircleOutlined style={{ color: DesignTokens.colors.info.main }} />
    }
  }

  const getTimeAgo = (timestamp: string): string => {
    const date = new Date(timestamp)
    const now = new Date()
    const diffMs = now.getTime() - date.getTime()
    const diffMins = Math.floor(diffMs / 60000)
    const diffHours = Math.floor(diffMs / 3600000)
    const diffDays = Math.floor(diffMs / 86400000)

    if (diffMins < 1) return 'Just now'
    if (diffMins < 60) return `${diffMins} minutes ago`
    if (diffHours < 24) return `${diffHours} hours ago`
    if (diffDays === 1) return 'Yesterday'
    if (diffDays < 7) return `${diffDays} days ago`
    return date.toLocaleDateString()
  }

  const dropdownContent = (
    <div className="notification-dropdown-content">
      <div className="notification-header">
        <Space style={{ justifyContent: 'space-between', width: '100%' }}>
          <Text strong>Notifications</Text>
          {unreadCount > 0 && (
            <Button
              type="link"
              size="small"
              icon={<CheckOutlined />}
              onClick={handleMarkAllAsRead}
              loading={markAllRead.isPending}
            >
              Mark all read
            </Button>
          )}
        </Space>
      </div>

      <div className="notification-list-container">
        {isLoading ? (
          <div style={{ textAlign: 'center', padding: '24px' }}>
            <Spin />
          </div>
        ) : notifications.length === 0 ? (
          <Empty
            image={Empty.PRESENTED_IMAGE_SIMPLE}
            description="No notifications"
            style={{ padding: '24px' }}
          />
        ) : (
          <List
            dataSource={notifications.slice(0, 5)}
            renderItem={(item) => (
              <List.Item
                className={`notification-item ${!item.read ? 'unread' : ''}`}
                onClick={() => handleNotificationClick(item)}
                actions={[
                  !item.read && (
                    <Button
                      key="mark-read"
                      type="text"
                      size="small"
                      icon={<CheckOutlined />}
                      onClick={(e) => handleMarkAsRead(item.id, e)}
                      title="Mark as read"
                    />
                  ),
                  <Button
                    key="delete"
                    type="text"
                    size="small"
                    danger
                    icon={<DeleteOutlined />}
                    onClick={(e) => handleDelete(item.id, e)}
                    title="Delete"
                  />,
                ].filter(Boolean)}
              >
                <List.Item.Meta
                  avatar={getNotificationIcon(item.type)}
                  title={
                    <div className="notification-title">
                      {item.title}
                      {!item.read && <span className="unread-dot" />}
                    </div>
                  }
                  description={
                    <div className="notification-description">
                      <Text type="secondary" style={{ fontSize: 12 }}>
                        {item.message}
                      </Text>
                      <Text type="secondary" style={{ fontSize: 11, marginTop: 4 }}>
                        {getTimeAgo((item as any).created_at || item.timestamp)}
                      </Text>
                    </div>
                  }
                />
              </List.Item>
            )}
          />
        )}
      </div>

      <div className="notification-footer">
        <Button
          type="link"
          block
          icon={<EyeOutlined />}
          onClick={() => {
            navigate('/notifications')
            setDropdownVisible(false)
          }}
        >
          View all notifications
        </Button>
      </div>
    </div>
  )

  return (
    <Dropdown
      dropdownRender={() => dropdownContent}
      trigger={['click']}
      placement="bottomRight"
      open={dropdownVisible}
      onOpenChange={setDropdownVisible}
    >
      <Badge count={unreadCount} showZero={false} className={className}>
        <BellOutlined className="notification-bell-icon" />
      </Badge>
    </Dropdown>
  )
}
