import { useState, useEffect, useCallback } from 'react'
import {
  Card,
  Row,
  Col,
  Avatar,
  Typography,
  Button,
  Statistic,
  Upload,
  Form,
  Input,
  Modal,
  message,
  Timeline,
  Tag,
  Divider,
  Space,
  Spin,
} from 'antd'
import {
  UserOutlined,
  EditOutlined,
  LockOutlined,
  MailOutlined,
  PhoneOutlined,
  CalendarOutlined,
  FolderOutlined,
  ExperimentOutlined,
  ClockCircleOutlined,
  CheckCircleOutlined,
  CloseCircleOutlined,
  UploadOutlined,
  WarningOutlined,
  InfoCircleOutlined,
} from '@ant-design/icons'
import { authStore } from '@/store/authStore'
import { PageHeader, StatisticCard } from '@/components/common'
import { DesignTokens } from '@/styles/design-tokens'
import { userService, type UserStats, type UserActivity } from '@/services/user.service'
import { logger } from '@/utils/logger'
import './ProfilePage.css'

const { Title, Text, Paragraph } = Typography

export const ProfilePage: React.FC = () => {
  const { user, refreshUser } = authStore()
  const [editModalVisible, setEditModalVisible] = useState(false)
  const [passwordModalVisible, setPasswordModalVisible] = useState(false)
  const [form] = Form.useForm()
  const [passwordForm] = Form.useForm()
  const [loading, setLoading] = useState(false)
  const [statsLoading, setStatsLoading] = useState(false)
  const [activitiesLoading, setActivitiesLoading] = useState(false)

  // User statistics
  const [stats, setStats] = useState<UserStats>({
    items: 0,
    samples: 0,
    tasks: 0,
    storageUsed: 0,
    storageTotal: 100,
  })

  // User activities
  const [activities, setActivities] = useState<UserActivity[]>([])

  // Fetch user stats
  const fetchUserStats = useCallback(async () => {
    setStatsLoading(true)
    try {
      const userStats = await userService.getUserStats()
      setStats(userStats)
    } catch (error) {
      logger.error('Failed to fetch user stats:', error)
    } finally {
      setStatsLoading(false)
    }
  }, [])

  // Fetch user activities
  const fetchUserActivities = useCallback(async () => {
    setActivitiesLoading(true)
    try {
      const userActivities = await userService.getUserActivity(10)
      setActivities(userActivities)
    } catch (error) {
      logger.error('Failed to fetch user activities:', error)
    } finally {
      setActivitiesLoading(false)
    }
  }, [])

  useEffect(() => {
    fetchUserStats()
    fetchUserActivities()
  }, [fetchUserStats, fetchUserActivities])

  const handleEditProfile = async (values: any) => {
    setLoading(true)
    try {
      await userService.updateProfile({
        username: values.username,
        email: values.email,
        phone: values.phone,
      })
      message.success('Profile updated successfully!')
      setEditModalVisible(false)
      form.resetFields()
      // Refresh user data in store
      refreshUser()
    } catch (error) {
      logger.error('Failed to update profile:', error)
      message.error('Failed to update profile')
    } finally {
      setLoading(false)
    }
  }

  const handleChangePassword = async (values: any) => {
    setLoading(true)
    try {
      await userService.changePassword({
        currentPassword: values.currentPassword,
        newPassword: values.newPassword,
      })
      message.success('Password changed successfully!')
      setPasswordModalVisible(false)
      passwordForm.resetFields()
    } catch (error) {
      logger.error('Failed to change password:', error)
      message.error('Failed to change password')
    } finally {
      setLoading(false)
    }
  }

  const handleAvatarUpload = async (file: File) => {
    try {
      await userService.uploadAvatar(file)
      message.success('Avatar updated successfully!')
      // Refresh user data to show new avatar
      refreshUser()
      return false // Prevent default upload behavior
    } catch (error) {
      logger.error('Failed to upload avatar:', error)
      message.error('Failed to upload avatar')
      return false
    }
  }

  const getActivityIcon = (type: string) => {
    switch (type) {
      case 'success':
        return <CheckCircleOutlined style={{ color: DesignTokens.colors.success.main }} />
      case 'error':
        return <CloseCircleOutlined style={{ color: DesignTokens.colors.error.main }} />
      case 'warning':
        return <WarningOutlined style={{ color: DesignTokens.colors.warning.main }} />
      case 'info':
        return <InfoCircleOutlined style={{ color: DesignTokens.colors.info.main }} />
      default:
        return <ClockCircleOutlined style={{ color: DesignTokens.colors.info.main }} />
    }
  }

  return (
    <div className="profile-page">
      <PageHeader title="个人资料" subtitle="管理您的账户信息和偏好设置" />

      <Row gutter={[24, 24]}>
        {/* Left Column: Profile Info */}
        <Col xs={24} lg={8}>
          <Card className="profile-card">
            <div className="profile-header">
              <div className="avatar-container">
                <Avatar size={120} icon={<UserOutlined />} className="profile-avatar" />
                <Upload accept="image/*" showUploadList={false} beforeUpload={handleAvatarUpload}>
                  <Button type="primary" size="small" icon={<UploadOutlined />} className="avatar-upload-btn">
                    更换头像
                  </Button>
                </Upload>
              </div>

              <Title level={3} className="profile-name">
                {user?.username}
              </Title>

              <Tag color={user?.role === 'admin' ? 'red' : 'blue'} className="role-tag">
                {user?.role === 'admin' ? '管理员' : '用户'}
              </Tag>
            </div>

            <Divider />

            <div className="profile-info">
              <div className="info-item">
                <MailOutlined className="info-icon" />
                <div>
                  <Text type="secondary" className="info-label">
                    邮箱
                  </Text>
                  <Text className="info-value">{user?.email || '未设置'}</Text>
                </div>
              </div>

              <div className="info-item">
                <PhoneOutlined className="info-icon" />
                <div>
                  <Text type="secondary" className="info-label">
                    电话
                  </Text>
                  <Text className="info-value">未设置</Text>
                </div>
              </div>

              <div className="info-item">
                <CalendarOutlined className="info-icon" />
                <div>
                  <Text type="secondary" className="info-label">
                    加入时间
                  </Text>
                  <Text className="info-value">2024-01-15</Text>
                </div>
              </div>
            </div>

            <Divider />

            <Space direction="vertical" style={{ width: '100%' }} size="middle">
              <Button type="primary" block icon={<EditOutlined />} onClick={() => setEditModalVisible(true)}>
                编辑资料
              </Button>
              <Button block icon={<LockOutlined />} onClick={() => setPasswordModalVisible(true)}>
                修改密码
              </Button>
            </Space>
          </Card>
        </Col>

        {/* Right Column: Statistics & Activity */}
        <Col xs={24} lg={16}>
          {/* Statistics */}
          <Spin spinning={statsLoading}>
            <Row gutter={[16, 16]}>
              <Col xs={24} sm={12} lg={6}>
                <StatisticCard
                  title="项目数"
                  value={stats.items}
                  icon={<FolderOutlined />}
                  color={DesignTokens.colors.primary.main}
                />
              </Col>
              <Col xs={24} sm={12} lg={6}>
                <StatisticCard
                  title="样本数"
                  value={stats.samples}
                  icon={<ExperimentOutlined />}
                  color={DesignTokens.colors.success.main}
                />
              </Col>
              <Col xs={24} sm={12} lg={6}>
                <StatisticCard
                  title="任务数"
                  value={stats.tasks}
                  icon={<ClockCircleOutlined />}
                  color={DesignTokens.colors.warning.main}
                />
              </Col>
              <Col xs={24} sm={12} lg={6}>
                <Card className="storage-card">
                  <Statistic
                    title="存储使用"
                    value={stats.storageUsed}
                    suffix={`/ ${stats.storageTotal} GB`}
                    precision={1}
                  />
                  <div
                    className="storage-bar"
                    style={{
                      width: '100%',
                      height: 8,
                      backgroundColor: DesignTokens.colors.gray[200],
                      borderRadius: 4,
                      overflow: 'hidden',
                      marginTop: 12,
                    }}
                  >
                    <div
                      className="storage-bar-fill"
                      style={{
                        width: `${(stats.storageUsed / stats.storageTotal) * 100}%`,
                        height: '100%',
                        backgroundColor:
                          stats.storageUsed / stats.storageTotal > 0.8
                            ? DesignTokens.colors.error.main
                            : DesignTokens.colors.primary.main,
                        transition: 'width 0.3s ease',
                      }}
                    />
                  </div>
                </Card>
              </Col>
            </Row>
          </Spin>

          {/* Activity Timeline */}
          <Card title="最近活动" className="activity-card" style={{ marginTop: 24 }}>
            {activitiesLoading ? (
              <div style={{ textAlign: 'center', padding: 24 }}>
                <Spin />
              </div>
            ) : activities.length === 0 ? (
              <div style={{ textAlign: 'center', padding: 24, color: DesignTokens.colors.gray[500] }}>
                No recent activities
              </div>
            ) : (
              <Timeline
                items={activities.map((activity) => ({
                  dot: getActivityIcon(activity.type),
                  children: (
                    <div className="activity-item">
                      <div className="activity-content">
                        <Text strong>{activity.title}</Text>
                        <Paragraph type="secondary" style={{ marginBottom: 0, marginTop: 4 }}>
                          {activity.description}
                        </Paragraph>
                      </div>
                      <Text type="secondary" className="activity-time">
                        {activity.time}
                      </Text>
                    </div>
                  ),
                }))}
              />
            )}
          </Card>
        </Col>
      </Row>

      {/* Edit Profile Modal */}
      <Modal
        title="编辑个人资料"
        open={editModalVisible}
        onCancel={() => {
          setEditModalVisible(false)
          form.resetFields()
        }}
        onOk={() => form.submit()}
        okText="保存"
        cancelText="取消"
        confirmLoading={loading}
      >
        <Form
          form={form}
          layout="vertical"
          onFinish={handleEditProfile}
          initialValues={{
            username: user?.username,
            email: user?.email,
          }}
        >
          <Form.Item label="用户名" name="username" rules={[{ required: true, message: '请输入用户名' }]}>
            <Input prefix={<UserOutlined />} placeholder="请输入用户名" />
          </Form.Item>

          <Form.Item
            label="邮箱"
            name="email"
            rules={[
              { required: true, message: '请输入邮箱' },
              { type: 'email', message: '请输入有效的邮箱地址' },
            ]}
          >
            <Input prefix={<MailOutlined />} placeholder="请输入邮箱" />
          </Form.Item>

          <Form.Item label="电话" name="phone">
            <Input prefix={<PhoneOutlined />} placeholder="请输入电话号码" />
          </Form.Item>
        </Form>
      </Modal>

      {/* Change Password Modal */}
      <Modal
        title="修改密码"
        open={passwordModalVisible}
        onCancel={() => {
          setPasswordModalVisible(false)
          passwordForm.resetFields()
        }}
        onOk={() => passwordForm.submit()}
        okText="确认修改"
        cancelText="取消"
        confirmLoading={loading}
      >
        <Form form={passwordForm} layout="vertical" onFinish={handleChangePassword}>
          <Form.Item label="当前密码" name="currentPassword" rules={[{ required: true, message: '请输入当前密码' }]}>
            <Input.Password prefix={<LockOutlined />} placeholder="请输入当前密码" />
          </Form.Item>

          <Form.Item
            label="新密码"
            name="newPassword"
            rules={[
              { required: true, message: '请输入新密码' },
              { min: 6, message: '密码至少6个字符' },
            ]}
          >
            <Input.Password prefix={<LockOutlined />} placeholder="请输入新密码" />
          </Form.Item>

          <Form.Item
            label="确认新密码"
            name="confirmPassword"
            dependencies={['newPassword']}
            rules={[
              { required: true, message: '请确认新密码' },
              ({ getFieldValue }) => ({
                validator(_, value) {
                  if (!value || getFieldValue('newPassword') === value) {
                    return Promise.resolve()
                  }
                  return Promise.reject(new Error('两次输入的密码不一致'))
                },
              }),
            ]}
          >
            <Input.Password prefix={<LockOutlined />} placeholder="请再次输入新密码" />
          </Form.Item>
        </Form>
      </Modal>
    </div>
  )
}
