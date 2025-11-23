import { useState, useEffect } from 'react'
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
import { DesignTokens } from '@/styles/design-tokens'
import './NotificationDropdown.css'

const { Text } = Typography

interface NotificationDropdownProps {
  className?: string
}

export const NotificationDropdown: React.FC<NotificationDropdownProps> = ({ className }) => {
  const navigate = useNavigate()
  const [notifications, setNotifications] = useState<Notification[]>([])
  const [loading, setLoading] = useState(false)
  const [unreadCount, setUnreadCount] = useState(0)
  const [dropdownVisible, setDropdownVisible] = useState(false)

  // Mock data - replace with actual API calls
  useEffect(() => {
    fetchNotifications()
  }, [])

  const fetchNotifications = async () => {
    setLoading(true)
    try {
      // TODO: Replace with actual API call
      // const response = await notificationService.getNotifications({ pageSize: 5 })
      // setNotifications(response.notifications)

      // Mock data
      const mockNotifications: Notification[] = [
        {
          id: '1',
          type: 'success',
          category: 'pipeline',
          title: 'Pipeline completed successfully',
          message: 'RNA-seq analysis pipeline has finished processing',
          timestamp: new Date(Date.now() - 2 * 60 * 60 * 1000).toISOString(),
          read: false,
          actionUrl: '/results/123',
          actionText: 'View Results',
        },
        {
          id: '2',
          type: 'error',
          category: 'task',
          title: 'Task failed',
          message: 'Alignment task encountered an error',
          timestamp: new Date(Date.now() - 5 * 60 * 60 * 1000).toISOString(),
          read: false,
          actionUrl: '/tasks/456',
          actionText: 'View Details',
        },
        {
          id: '3',
          type: 'info',
          category: 'system',
          title: 'System update',
          message: 'New features are now available',
          timestamp: new Date(Date.now() - 24 * 60 * 60 * 1000).toISOString(),
          read: true,
        },
        {
          id: '4',
          type: 'warning',
          category: 'project',
          title: 'Storage limit warning',
          message: 'Project storage is 85% full',
          timestamp: new Date(Date.now() - 2 * 24 * 60 * 60 * 1000).toISOString(),
          read: true,
          actionUrl: '/items',
          actionText: 'Manage Storage',
        },
      ]
      setNotifications(mockNotifications)
      setUnreadCount(mockNotifications.filter((n) => !n.read).length)
    } catch (error) {
      console.error('Failed to fetch notifications:', error)
    } finally {
      setLoading(false)
    }
  }

  const handleMarkAsRead = async (id: string, event: React.MouseEvent) => {
    event.stopPropagation()
    try {
      // TODO: Call API to mark as read
      // await notificationService.markAsRead(id)
      setNotifications(
        notifications.map((n) => (n.id === id ? { ...n, read: true } : n))
      )
      setUnreadCount((prev) => Math.max(0, prev - 1))
    } catch (error) {
      console.error('Failed to mark notification as read:', error)
    }
  }

  const handleMarkAllAsRead = async () => {
    try {
      // TODO: Call API to mark all as read
      // await notificationService.markAllAsRead()
      setNotifications(notifications.map((n) => ({ ...n, read: true })))
      setUnreadCount(0)
    } catch (error) {
      console.error('Failed to mark all as read:', error)
    }
  }

  const handleDelete = async (id: string, event: React.MouseEvent) => {
    event.stopPropagation()
    try {
      // TODO: Call API to delete notification
      // await notificationService.deleteNotification(id)
      const notification = notifications.find((n) => n.id === id)
      setNotifications(notifications.filter((n) => n.id !== id))
      if (notification && !notification.read) {
        setUnreadCount((prev) => Math.max(0, prev - 1))
      }
    } catch (error) {
      console.error('Failed to delete notification:', error)
    }
  }

  const handleNotificationClick = (notification: Notification) => {
    if (!notification.read) {
      handleMarkAsRead(notification.id, {} as React.MouseEvent)
    }
    if (notification.actionUrl) {
      navigate(notification.actionUrl)
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
            >
              Mark all read
            </Button>
          )}
        </Space>
      </div>

      <div className="notification-list-container">
        {loading ? (
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
                      type="text"
                      size="small"
                      icon={<CheckOutlined />}
                      onClick={(e) => handleMarkAsRead(item.id, e)}
                      title="Mark as read"
                    />
                  ),
                  <Button
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
                        {getTimeAgo(item.timestamp)}
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
