/**
 * Main Layout - For authenticated pages
 */
import React, { useState } from 'react'
import { Outlet, useNavigate, useLocation } from 'react-router-dom'
import { Layout, Menu, Typography, Dropdown, Avatar, Space } from 'antd'
import {
  DashboardOutlined,
  FolderOutlined,
  ExperimentOutlined,
  BarChartOutlined,
  SettingOutlined,
  UserOutlined,
  LogoutOutlined,
  MenuFoldOutlined,
  MenuUnfoldOutlined,
  FileOutlined,
  ClockCircleOutlined,
  ThunderboltOutlined,
  RobotOutlined,
  FundOutlined,
  BookOutlined,
} from '@ant-design/icons'
import { authStore } from '@/store/authStore'
import { ThemeToggle } from '@/components/common'
import { NotificationDropdown } from '@/components/notifications/NotificationDropdown'
import styles from './MainLayout.module.css'

const { Header, Sider, Content } = Layout
const { Text } = Typography

export const MainLayout: React.FC = () => {
  const navigate = useNavigate()
  const location = useLocation()
  const { user, logout } = authStore()
  const [collapsed, setCollapsed] = useState(false)

  const menuItems = [
    {
      key: '/dashboard',
      icon: <DashboardOutlined />,
      label: 'Dashboard',
    },
    {
      key: '/ai',
      icon: <RobotOutlined />,
      label: 'AI Intelligence',
    },
    {
      key: '/items',
      icon: <FolderOutlined />,
      label: 'Projects',
    },
    {
      key: '/samples',
      icon: <ExperimentOutlined />,
      label: 'Samples',
    },
    {
      key: '/files',
      icon: <FileOutlined />,
      label: 'Files',
    },
    {
      key: '/pipelines',
      icon: <ThunderboltOutlined />,
      label: 'Pipelines',
    },
    {
      key: '/tasks',
      icon: <ClockCircleOutlined />,
      label: 'Tasks',
    },
    {
      key: '/results',
      icon: <BarChartOutlined />,
      label: 'Results',
    },
    {
      key: '/analytics',
      icon: <FundOutlined />,
      label: 'Analytics',
    },
    {
      key: '/knowledge',
      icon: <BookOutlined />,
      label: 'Knowledge',
    },
    ...(user?.role === 'admin'
      ? [
          {
            key: '/admin',
            icon: <SettingOutlined />,
            label: 'Admin',
          },
        ]
      : []),
  ]

  const userMenuItems = [
    {
      key: 'profile',
      icon: <UserOutlined />,
      label: 'Profile',
      onClick: () => navigate('/profile'),
    },
    {
      key: 'settings',
      icon: <SettingOutlined />,
      label: 'Settings',
      onClick: () => navigate('/settings'),
    },
    {
      type: 'divider' as const,
    },
    {
      key: 'logout',
      icon: <LogoutOutlined />,
      label: 'Logout',
      onClick: () => {
        logout()
        navigate('/login')
      },
    },
  ]

  const handleMenuClick = ({ key }: { key: string }) => {
    navigate(key)
  }

  return (
    <Layout className={styles.layout}>
      <Sider trigger={null} collapsible collapsed={collapsed} className={styles.sider} width={240}>
        <div className={styles.logo}>
          <ExperimentOutlined style={{ fontSize: 24, color: '#2196F3' }} />
          {!collapsed && (
            <Text strong style={{ color: '#fff', marginLeft: 12, fontSize: 16 }}>
              NGSmodule
            </Text>
          )}
        </div>

        <Menu
          theme="dark"
          mode="inline"
          selectedKeys={[location.pathname]}
          items={menuItems}
          onClick={handleMenuClick}
        />
      </Sider>

      <Layout>
        <Header className={styles.header}>
          <div className={styles.headerLeft}>
            {React.createElement(collapsed ? MenuUnfoldOutlined : MenuFoldOutlined, {
              className: styles.trigger,
              onClick: () => setCollapsed(!collapsed),
            })}
          </div>

          <div className={styles.headerRight}>
            <Space size="large">
              <ThemeToggle mode="icon" size="default" />

              <NotificationDropdown className={styles.headerIcon} />

              <Dropdown menu={{ items: userMenuItems }} placement="bottomRight">
                <Space className={styles.userInfo}>
                  <Avatar icon={<UserOutlined />} style={{ backgroundColor: '#2196F3' }} />
                  <div className={styles.userDetails}>
                    <Text strong>{user?.username}</Text>
                    <Text type="secondary" style={{ fontSize: 12 }}>
                      {user?.role === 'admin' ? 'Administrator' : 'User'}
                    </Text>
                  </div>
                </Space>
              </Dropdown>
            </Space>
          </div>
        </Header>

        <Content className={styles.content}>
          <div className={styles.contentInner}>
            <Outlet />
          </div>
        </Content>
      </Layout>
    </Layout>
  )
}

export default MainLayout
